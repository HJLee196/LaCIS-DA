rm(list=ls())

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

multiknk = function(W_imp, alpha=0.1){
  m = dim(W_imp)[1] - 1
  
  W_max_ind = apply(W_imp,2,which.max)
  W_max = sapply(1:length(W_max_ind), function(x) W_imp[W_max_ind[x],x])
  W_med = sapply(1:length(W_max_ind), function(x) median(W_imp[-W_max_ind[x],x]))
  
  tau_stab <- (W_max - W_med)
  kap = (W_max_ind==1)
  
  ts = sort(c(0, abs(tau_stab)))
  ratio = sapply(ts, function(t) (1/m + 1/m*sum((tau_stab >= t)*(1-kap)))/max(1, 
                                                                              sum((tau_stab >= t)*(kap))))
  ok = which(ratio <= alpha)
  threshold <- ifelse(length(which(ts[ok]>0)) > 0, ts[ok][which(ts[ok]>0)[1]], Inf)
  
  selected <- which(tau_stab*kap >= threshold)
  return(list(rej = selected, cut = threshold))
}

# PC: Number of principal components that will be used for imputation and covariance estimation

seed_num = 17
set.seed(seed_num)
down = NULL
gene.index = NULL
m_kos = 10 # the number of the sets of knockoffs
llam = 0.1 
nSim = 1
#testUse = "MAST"
celltype = "microglia" #"astrocyte"
PC = 11 #43 #60 #115
ko_method = "multidecomp"
ko_test = "LR" #"mast"
lfc = c(0.0, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.15, 0.20, 0.25, 0.5)
#"wilcox"
#"LR"
#"MAST"
#"LCD"
sign_strength=3
n_signal = 50
hyperparam_used_list = 
  list(seed_num = seed_num,
       down = down,
       PC = PC,
       llam = llam,
       nSim = nSim)

#################################################################################
#  Load External Data
#################################################################################
# WARNING: The wrapper function 'Seurat_repeat' is written for this specific dataset. Changes will be required if 
# using a different dataset.

datafilename <- paste0("./HuVascAD_", celltype,".rds")
HuVascAD <- readRDS(file = datafilename)
options(future.globals.maxSize= +Inf)

if(is.null(down)){
  test_data <- HuVascAD # use full dataset
} else {
  healthy_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "Control"] 
  ad_M <- names(HuVascAD$treat)[HuVascAD$treat %in% "AD"] 
  
  # the function subset() might have a set.seed function somewhere, 
  # returning non-random results after the first run.
  set.seed(seed_num)
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
xp_data <- GetAssayData(object = xp_data.sub, layer = "count")
# rotate to match downstream analyses.
xp_data.matrix.exp <- t(as.matrix(xp_data) > 0) # indicator matrix for expressed genes, 
xp_data.matrix.exp.count <- apply(xp_data.matrix.exp,2,sum)

# xp_data.matrix <- (as.matrix(xp_data)) # do not transpose here. Keep the genes as rows.


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
# rm(test_data) # release memory space

xp_data <- GetAssayData(object = xp_data.sub, layer = "data") # Normalized data matrix


xp_data.matrix <- (as.matrix(xp_data))
xp_data.matrix <- t(xp_data.matrix)

print(dim(xp_data.matrix))
print(dim(xp_data.matrix.exp))

xp_data.sub$age2 = xp_data.sub$age^2
xp_data.sub$CDR = CDR

n <- dim(xp_data.matrix)[1]
#for(i in 1:length(err)){
#  xp_data.matrix.imp[,i] <-  xp_data.matrix.imp[,i] + rnorm(n,0,sqrt(err[i]))*(1-xp_data.matrix.exp[,i])
#}

end.time = Sys.time()

### Add covariates, we perform the test conditional on knockoff-selected variables

ko_result_filename <- paste0("HuVascAD_cell", celltype,
                             "_method", ko_method,
                             "_testUse", ko_test,
                             "_downNULL_gene23537_PC", PC,
                             "_10llam1_seed2017_maxIterImp100.rda")

load(file = ko_result_filename)

W_imp = FinalResult[[1]]$W_imp_b

m = 10 # since we will fix M at 10.
q = 0.1

# Average knockoff values
K_imp <- matrix(NA, m, dim(W_imp)[2])
for(w_it in 1:m){
  K_imp[w_it,] <- W_imp[1,] - W_imp[(w_it+1),]
}

ko_q<- multiknk(W_imp,q)
ko_selected_ind = multiknk(W_imp,q)$rej
ko_selected <- FinalResult[[1]]$feature.names.new[ko_selected_ind]
ko_selected

for (l in 1:length(ko_selected)) {
  varname <- paste0("ko_", l)
  xp_data.sub[[varname]] <- xp_data.matrix[, ko_selected_ind[l]]
}

remove(FinalResult)

### BH-selected genes

bh_resultname = paste0("HuVascAD_cell", celltype,
                       "_testUseBH_downNULL_gene23537_PC", PC,
                       "_10llam1_seed17_maxIterImp100.rda")
  
load(file = bh_resultname)

mast_p_value <- FinalResult[[1]]$selection_result_mast$p_val_adj

bh_selected_ind <- order(mast_p_value)[1:length(ko_selected)]

bh_selected<- FinalResult[[1]]$feature.names.new[bh_selected_ind]
bh_selected

for (l in 1:length(bh_selected)) {
  varname <- paste0("bh_", l)
  xp_data.sub[[varname]] <- xp_data.matrix[, bh_selected_ind[l]]
}

### random selected index 
random_ind <- sample.int(length(FinalResult[[1]]$feature.names.new), length(ko_selected))
random_selected <- FinalResult[[1]]$feature.names.new[random_ind]
random_selected
for (l in 1:length(random_selected)) {
  varname <- paste0("random_", l)
  xp_data.sub[[varname]] <- xp_data.matrix[, random_ind[l]]
}

rm(FinalResult)

FinalResult = list()

start_time_nsim_set = NULL
end_time_nsim_set = NULL

for(nsim in 1:nSim){
  # Computing Time
  start_time_nsim = Sys.time()
  start_time_nsim_set = c(start_time_nsim_set, start_time_nsim)
  
  set.seed(seed = seed_num + nsim)
  
  # Set labels to be observed ones
  Idents(xp_data.sub) <- 'treat'
  
  ####
  # 4-4-1. Variable Selection 
  ####
  
  ### conditional on ko-selected genes
  ### LR statistics
  q_value_mast_cond <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  
  xp_data.sub$batch <- as.factor(xp_data.sub$batch)
  xp_data.sub$gender<- as.factor(xp_data.sub$gender)
  covariates_list <- c(c("gender","age","age2","batch","CDR"), paste0("ko_", c(1:length(ko_selected))))
  result.cond.mast <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                  min.pct = 0,logfc.threshold=0,verbose = FALSE, test.use ="LR",
                                  latent.vars = covariates_list, mean.fxn = mean_fxn_v4)
  
  result.cond.mast.reordered <- result.cond.mast[match(feature.names.new, rownames(result.cond.mast)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.cond.mast.reordered$avg_log2FC) >= thres)
    tmp_p <- result.cond.mast.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_mast_cond[i, tmp_index] <- tmp_q
    q_value_mast_cond[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_mast_cond) <- lfc
  
  ### Conditional on BH-selected genes
  ### LR-statistics 
  
  q_value_mast_bh <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  
  covariates_list <- c(c("gender","age","age2","batch","CDR"), paste0("bh_",c(1:length(bh_selected))))
  result.bh.mast <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                min.pct = 0,logfc.threshold=0,verbose = TRUE, test.use ="LR",
                                latent.vars = covariates_list, mean.fxn = mean_fxn_v4)
  
  result.bh.mast.reordered <- result.bh.mast[match(feature.names.new, rownames(result.bh.mast)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.bh.mast.reordered$avg_log2FC) >= thres)
    tmp_p <- result.bh.mast.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_mast_bh[i, tmp_index] <- tmp_q
    q_value_mast_bh[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_mast_bh) <- lfc
  
  ### Conditional on random-selected genes
  ### LR-statistics 
  
  q_value_mast_rand <- matrix(NA,nrow = length(lfc), ncol = length(feature.names.new))
  
  covariates_list <- c(c("gender","age","age2","batch","CDR"), paste0("random_", c(1:length(random_selected))))
  result.rand.mast <- FindMarkers(xp_data.sub, ident.1 = 'Control', ident.2 = 'AD', slot = "data",
                                  min.pct = 0,logfc.threshold=0,verbose = TRUE, test.use ="LR",
                                  latent.vars = covariates_list, mean.fxn = mean_fxn_v4)
  
  result.rand.mast.reordered <- result.rand.mast[match(feature.names.new, rownames(result.rand.mast)),]
  
  for(i in 1:length(lfc)) {
    thres = lfc[i]
    tmp_index <- which(abs(result.rand.mast.reordered$avg_log2FC) >= thres)
    tmp_p <- result.rand.mast.reordered$p_val[tmp_index]
    tmp_q <- p.adjust(tmp_p, method="BH")
    q_value_mast_rand[i, tmp_index] <- tmp_q
    q_value_mast_rand[i, -tmp_index] <- 2 # set the genes with lfc below the threshold with q-values 2. 
  }
  rownames(q_value_mast_rand) <- lfc
  
  cat("# Iter: ", nsim, "\n")
  
  # Computing Time
  end_time_nsim = Sys.time()
  end_time_nsim_set = c(end_time_nsim_set, end_time_nsim)
  
  nsim_total_seconds = as.numeric(difftime(end_time_nsim , start_time_nsim, units = "secs"))
  nsim_minutes <- floor(nsim_total_seconds / 60)
  nsim_seconds <- round(nsim_total_seconds %% 60)
  
  print(sprintf("Spent Computation Time (each sim): %02d:%02d", nsim_minutes, nsim_seconds))
  
  FinalResult[[nsim]] <- (list(selection_result_mast_cond = result.cond.mast.reordered, 
                               q_value_mast_cond = q_value_mast_cond, 
                               ko_selected_ind = ko_selected_ind,
                               selection_result_mast_bh = result.bh.mast.reordered, 
                               q_value_mast_bh = q_value_mast_bh, 
                               bh_selected_ind = bh_selected_ind,
                               selection_result_mast_rand = result.rand.mast.reordered, 
                               q_value_mast_rand = q_value_mast_rand, 
                               random_ind = random_ind,
                               #A=A,
                               #B=B,
                               feature.names.new = feature.names.new,
                               hyperparam_used_list = hyperparam_used_list))
}
#################################################################################
#  Running the function nSim times:
#################################################################################


#################################################################################
#  Save the result to a folder on the cluster:
#################################################################################

outputname =
  bh_resultname %>%
  str_replace(pattern = "\\.rda", replacement = str_c("_conditional_", ko_method, "_", ko_test, "\\.rda")) 

save(FinalResult,file = outputname)

