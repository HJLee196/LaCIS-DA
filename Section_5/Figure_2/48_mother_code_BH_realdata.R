rm(list=ls())

###
# 1. Preparation
###

###
# 1-1. Load Necessary R Packages
###

library(Seurat)
library(knockoff)
library(MAST)
library(Matrix)
library(softImpute)
library(tidyverse)

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

# Enforce Seurat v4-style logFC in Seurat v5: use a fixed pseudocount (+1)
mean_fxn_v4 <- function(x) {
  # x: data layer values (log1p scale), provided as a matrix [genes x cells]
  #     or possibly as a vector
  # Steps:
  # 1) Convert back to the original scale using expm1()
  # 2) Compute the mean expression across cells
  # 3) Add a fixed pseudocount (+1)
  # 4) Take log2 to obtain log fold change
  log(rowMeans(expm1(x)) + 1, base = 2)
}

####
# 1-2. Read the dataset
#### 
# "path/to/your/rds_folder"
file_path = "path/to/your/rds_folder"
cell_name = 'astrocyte'
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

seed_num = 17

set.seed(seed_num)
down = NULL
gene.index = 1:max_gene
n_signal = 50
PC = 43
llam = 0.1 
nSim = 1
m_kos = 10 # the number of the sets of knockoffs
max_iter_imp = 100
testUse = "BH"
lfc = c(0.0, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.20, 0.25, 0.5)

hyperparam_used_list = 
  list(seed_num = seed_num,
       down = down,
       max_gene.index = max(gene.index),
       # sign_strength = sign_strength,
       n_signal = n_signal,
       PC = PC,
       llam = llam,
       nSim = nSim,
       max_iter_imp = max_iter_imp)

# if ((LCD_crit_use == "lambda.1se")&(testUse == "LCD"))
# {
#   testUseFull = str_c(testUse, "1")
# } else {
  testUseFull = testUse
# }

file_name =
  str_c("HuVascAD",
        "_cell", cell_name,
        # "_method", method_name,
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

####
# 1-4. Downsample the data
####

if(is.null(down)){
  test_data <- HuVascAD # use full dataset
} else {
  
  healthy_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "Control"] 
  ad_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "AD"] 
  
  # set.seed(seed_num) # To reproduce the same results, commented
  down_ind <- c(sample(healthy_M,down),
                sample(ad_M,down))
  which(names(HuVascAD$treat) %in% down_ind)
  
  test_data <- subset(HuVascAD, 
                      cells = which(names(HuVascAD$treat) %in% down_ind))
}
rm(HuVascAD) # rm dataset to save memory

# Take a subset of genes
if(length(gene.index) > 0){
  feature.names <- rownames(test_data)[gene.index]
} else {
  feature.names <- rownames(test_data)
}

xp_data.sub <- subset(x = test_data, features = feature.names)
DefaultAssay(xp_data.sub) <- "RNA"

####
# 1-5. Create the X variables
####

#### Calculate CDR ####
xp_data <- GetAssayData(object = xp_data.sub, layer = "count")

xp_data.matrix.exp <- t(as.matrix(xp_data) > 0) # indicator matrix for expressed genes, 
xp_data.matrix.exp.count <- apply(xp_data.matrix.exp,2,sum)

# Prespecify how many variables are fixed during the imputation
k_num = 7

# exclude empty variables
feature.exp.ind <- which(xp_data.matrix.exp.count > (PC+k_num+1)) 
feature.names.new <- feature.names[feature.exp.ind]  # length(feature.names.new)

xp_data.matrix.exp <- xp_data.matrix.exp[,feature.exp.ind]
xp_data.matrix.exp.count <- xp_data.matrix.exp.count[feature.exp.ind]
CDR = rowMeans(xp_data.matrix.exp) # expressed genes in each CELL


# refresh subset
xp_data.sub <- subset(x = test_data, features = feature.names.new)
xp_data <- GetAssayData(object = xp_data.sub, layer = "data") # Normalized data matrix


xp_data.matrix <- (as.matrix(xp_data))
xp_data.matrix <- t(xp_data.matrix)

print(dim(xp_data.matrix))
print(dim(xp_data.matrix.exp))

####
# 2. Imputation for the factor variables
####

# Center target matrix, but only center the expressed parts.
xp_data.matrix.fix <- (colSums(xp_data.matrix, na.rm = T))/xp_data.matrix.exp.count
xp_data.matrix <- sweep(xp_data.matrix,2,xp_data.matrix.fix,FUN = "-")*(xp_data.matrix.exp)

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
colnames(batch) <- c("batch2", "batch3", "batch4")

xp_data.sub$age2 = xp_data.sub$age^2
xp_data.sub$CDR = CDR
xp_data.sub$batch2 <- batch[,1]
xp_data.sub$batch3 <- batch[,2]
xp_data.sub$batch4 <- batch[,3]

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

r = ncol(xp_data.matrix.imp_l[[2]])

A_sub <- xp_data.matrix.imp_l[[2]] [, (k_num+1):r]
B_sub <- t(xp_data.matrix.imp_l[[3]])[, (k_num+1):r]
B_rotated <- varimax(B_sub)

B_varimax <- B_rotated$loadings
A_varimax <- A_sub %*% B_rotated$rotmat

xp_data.sub$age2 = xp_data.sub$age^2

n <- dim(xp_data.matrix.imp)[1]
#for(i in 1:length(err)){
#  xp_data.matrix.imp[,i] <-  xp_data.matrix.imp[,i] + rnorm(n,0,sqrt(err[i]))*(1-xp_data.matrix.exp[,i])
#}

for (l in 1:PC) {
  varname <- paste0("Factor_", l)
  xp_data.sub[[varname]] <- A_varimax[,l]
}

end.time = Sys.time()

####
# 3. Run simulations
####

if ("test_data" %in% ls()){
  rm(test_data)
}

FinalResult = list()

start_time_nsim_set = NULL
end_time_nsim_set = NULL

for(nsim in 1:nSim){
  # Computing Time
  start_time_nsim = Sys.time()
  start_time_nsim_set = c(start_time_nsim_set, start_time_nsim)
  
  # As Yixia mentioned, some function might have set.seed in it.
  set.seed(seed = seed_num + nsim)
  
  # Set labels to be observed ones
  Idents(xp_data.sub) <- 'treat'
  
  ####
  # 4-4-1. Variable Selection 
  ####
  
  ### MAST statistic
  q_value_mast <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  
  # shouldn't matter which slot I change as long as I use it.
  # xp_data.sub <- SetAssayData(object = xp_data.sub, slot = "count", new.data = t(xp_data.matrix))
  xp_data.sub$batch <- as.factor(xp_data.sub$batch)
  xp_data.sub$gender<- as.factor(xp_data.sub$gender)
  result.sub.mast <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                            min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
                            latent.vars = c("gender","age","age2","batch","CDR"), mean.fxn = mean_fxn_v4)
  # result.sub_7covariate <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
  #                           min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
  #                           latent.vars = c("gender","age","age2","batch2","batch3","batch4", "CDR"))
  
  result.sub.mast.reordered <- result.sub.mast[match(feature.names.new, rownames(result.sub.mast)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.sub.mast.reordered$avg_log2FC) >= thres)
    tmp_p <- result.sub.mast.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_mast[i, tmp_index] <- tmp_q
    q_value_mast[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_mast) <- lfc
  
  ### Adjust with respect to covariates
  
  q_value_mast_adj <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  #p_value_b <- matrix(NA,nrow = 2, ncol = length(feature.names.new))
  
  # Calculate importance statistic for each knockoffs
  
  # use observed values to find p-values
  # shouldn't matter which slot I change as long as I use it.
  #xp_data.sub <- SetAssayData(object = xp_data.sub, slot = "count", new.data = t(xp_data.matrix))
  #xp_data.sub$batch <- as.factor(xp_data.sub$batch)
  #xp_data.sub$gender<- as.factor(xp_data.sub$gender)
  latent.vars.list <- c(c("gender","age","age2","batch","CDR"), paste0("Factor_", c(1:PC)))
  result.sub.mast.adj <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                 min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
                                 latent.vars = latent.vars.list, mean.fxn = mean_fxn_v4)
  # result.sub_7covariate <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
  #                           min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
  #                           latent.vars = c("gender","age","age2","batch2","batch3","batch4", "CDR"))
  
  result.sub.mast.adj.reordered <- result.sub.mast.adj[match(feature.names.new, rownames(result.sub.mast.adj)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.sub.mast.adj.reordered$avg_log2FC) >= thres)
    tmp_p <- result.sub.mast.adj.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_mast_adj[i, tmp_index] <- tmp_q
    q_value_mast_adj[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_mast_adj) <- lfc
  
  ### LRT statistic
  q_value_lr <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  
  #xp_data.sub <- SetAssayData(object = xp_data.sub, slot = "count", new.data = t(xp_data.matrix))
  #xp_data.sub$batch <- as.factor(xp_data.sub$batch)
  #xp_data.sub$gender<- as.factor(xp_data.sub$gender)
  result.sub.lr <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                 min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="LR",
                                 latent.vars = c("gender","age","age2","batch","CDR"), mean.fxn = mean_fxn_v4)
  # result.sub_7covariate <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
  #                           min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
  #                           latent.vars = c("gender","age","age2","batch2","batch3","batch4", "CDR"))
  
  result.sub.lr.reordered <- result.sub.lr[match(feature.names.new, rownames(result.sub.lr)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.sub.lr.reordered$avg_log2FC) >= thres)
    tmp_p <- result.sub.lr.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_lr[i, tmp_index] <- tmp_q
    q_value_lr[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_lr) <- lfc
  
  ### Adjust with respect to covariates
  
  q_value_lr_adj <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
 
  #xp_data.sub <- SetAssayData(object = xp_data.sub, slot = "count", new.data = t(xp_data.matrix))
  #xp_data.sub$batch <- as.factor(xp_data.sub$batch)
  #xp_data.sub$gender<- as.factor(xp_data.sub$gender)
  latent.vars.list <- c(c("gender","age","age2","batch","CDR"), paste0("Factor_", c(1:PC)))
  result.sub.lr.adj <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                     min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="LR",
                                     latent.vars = latent.vars.list, mean.fxn = mean_fxn_v4)
  # result.sub_7covariate <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
  #                           min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
  #                           latent.vars = c("gender","age","age2","batch2","batch3","batch4", "CDR"))
  
  result.sub.lr.adj.reordered <- result.sub.lr.adj[match(feature.names.new, rownames(result.sub.lr.adj)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.sub.lr.adj.reordered$avg_log2FC) >= thres)
    tmp_p <- result.sub.lr.adj.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_lr_adj[i, tmp_index] <- tmp_q
    q_value_lr_adj[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_lr_adj) <- lfc
  
  ### WRT statistic
  q_value_wrt <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  
  #xp_data.sub <- SetAssayData(object = xp_data.sub, slot = "count", new.data = t(xp_data.matrix))
  #xp_data.sub$batch <- as.factor(xp_data.sub$batch)
  #xp_data.sub$gender<- as.factor(xp_data.sub$gender)
  result.sub.wrt <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                               min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="wilcox",
                               latent.vars = NULL, mean.fxn = mean_fxn_v4)
  # result.sub_7covariate <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
  #                           min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
  #                           latent.vars = c("gender","age","age2","batch2","batch3","batch4", "CDR"))
  
  result.sub.wrt.reordered <- result.sub.wrt[match(feature.names.new, rownames(result.sub.wrt)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.sub.wrt.reordered$avg_log2FC) >= thres)
    tmp_p <- result.sub.wrt.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_wrt[i, tmp_index] <- tmp_q
    q_value_wrt[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_wrt) <- lfc
  
  ### Adjust with respect to covariates
  
  q_value_wrt_adj <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  
  #xp_data.sub <- SetAssayData(object = xp_data.sub, slot = "count", new.data = t(xp_data.matrix))
  #xp_data.sub$batch <- as.factor(xp_data.sub$batch)
  #xp_data.sub$gender<- as.factor(xp_data.sub$gender)
  latent.vars.list <- c(c("gender","age","age2","batch","CDR"), paste0("Factor_", c(1:PC)))
  result.sub.wrt.adj <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                   min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="wilcox",
                                   latent.vars = NULL, mean.fxn = mean_fxn_v4)
  # result.sub_7covariate <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
  #                           min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="MAST",
  #                           latent.vars = c("gender","age","age2","batch2","batch3","batch4", "CDR"))
  
  result.sub.wrt.adj.reordered <- result.sub.wrt.adj[match(feature.names.new, rownames(result.sub.wrt.adj)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.sub.wrt.adj.reordered$avg_log2FC) >= thres)
    tmp_p <- result.sub.wrt.adj.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_wrt_adj[i, tmp_index] <- tmp_q
    q_value_wrt_adj[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_wrt_adj) <- lfc
  
  cat("# Iter: ", nsim, "\n")
  
  # Computing Time
  end_time_nsim = Sys.time()
  end_time_nsim_set = c(end_time_nsim_set, end_time_nsim)
  
  nsim_total_seconds = as.numeric(difftime(end_time_nsim , start_time_nsim, units = "secs"))
  nsim_minutes <- floor(nsim_total_seconds / 60)
  nsim_seconds <- round(nsim_total_seconds %% 60)
  
  print(sprintf("Spent Computation Time (each sim): %02d:%02d", nsim_minutes, nsim_seconds))
  
  FinalResult[[nsim]] <- (list(selection_result_mast = result.sub.mast.reordered, 
                               q_value_mast = q_value_mast, 
                               selection_result_mast_adj = result.sub.mast.adj.reordered, 
                               q_value_mast_adj = q_value_mast_adj, 
                               selection_result_lr = result.sub.lr.reordered, 
                               q_value_lr = q_value_lr, 
                               selection_result_lr_adj = result.sub.lr.adj.reordered, 
                               q_value_lr_adj = q_value_lr_adj, 
                               selection_result_wrt = result.sub.wrt.reordered, 
                               q_value_wrt = q_value_wrt, 
                               selection_result_wrt_adj = result.sub.wrt.adj.reordered, 
                               q_value_wrt_adj = q_value_wrt_adj, 
                               feature.names.new = feature.names.new,
                               hyperparam_used_list = hyperparam_used_list))
}

####
# 5. Save the results
####

# outputname = paste0("HuVascAD_downNULL_geneNULL_BH_PC",PC, "_","hippo_onlyFALSE_",cell_name,"_.rda" )
# save(FinalResult,file = outputname)
save(FinalResult,
     file = file_name)

