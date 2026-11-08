library("stats")
#setwd("")

### Load the imputation/factor analysis results
### This is obtained by running the script "HuVascAD_downNULL_geneNULL_scsoftImpute_PC43_ll01_hippo_onlyFALSE_astrocyte_factors.R"
load("./factors/HuVascAD_downNULL_geneNULL_scsoftImpute_PC43_hippo_onlyFALSE_astrocyte_factors.rda")
A <- FinalResult[[1]]
B <- FinalResult[[2]]
n = nrow(A)
p = nrow(B)
r = ncol(B)
k = 7 # number of covariates used
PC = 43 # number of factors for astrocyte data

A_var <- apply(A, 2, var)
A_var
B_var <- apply(B, 2, var)
B_var

### Remove the covariates, extract the remaining latent factors
A_sub <- A[, (k+1):r]
B_sub <- B[, (k+1):r]

### Apply the varimax rotation
B_rotated <- varimax(B_sub)
B_varimax <- B_rotated$loadings
A_varimax <- A_sub %*% B_rotated$rotmat

colnames(B_varimax) <- paste0("Factor ", c(1:(r-k)))

###################################################################
###################################################################
###################################################################
### Generate histograms of factor loadings
library(ggplot2)
library(tidyr)
library(dplyr)

# Reshape the dataset to long format
plot_B <- matrix(B_varimax, nrow =p, ncol=PC)
plot_B <- as.data.frame(plot_B)
colnames(plot_B) <- paste0("Factor ", c(1:(r-k)))
plot_B <- plot_B %>% select(1:16)        # select only first 16 variables
loadings_long <- plot_B %>%
  pivot_longer(cols = everything(),  # Reshape all columns
               names_to = "variable",  # New column for variable names
               values_to = "value")    # New column for values

loadings_long$variable <- factor(loadings_long$variable, levels = paste0("Factor ", c(1:16)))

# Create a histogram for factor loadings
loading_plot<- ggplot(loadings_long, aes(x = value)) +
  geom_histogram(binwidth = 0.02, fill = "skyblue", color = NA) +
  facet_wrap(~ variable, scales = "free") +  # Create a panel for each variable
  labs(title = "Histograms of factor loadings",
       x = "Value", y = "Count") +
  theme(
    plot.title = element_text(size = 14),
    axis.title.x = element_text(size = 11),
    axis.title.y = element_text(size = 11)
  )

loading_plot

jpeg(filename = paste0("factor_loadings.jpeg"), width = 7, height = 7, units = "in",
     res = 600)
print(loading_plot)
dev.off()

###################################################################
###################################################################
###################################################################
### Code for exploratory analysis ###
### Generat UMAPs to understand meanings of each factor ###

suppressPackageStartupMessages({
  library(Seurat)
  #library(knockoff)
  #library(MAST)
  library(Matrix)
  library(softImpute)
  library(tidyverse)
  library(RColorBrewer)
})
#################################################################################
#  Load External Data
#################################################################################
# WARNING: The wrapper function 'Seurat_repeat' is written for this specific dataset. Changes will be required if 
# using a different dataset.
#setwd("")
#HuVascAD <- readRDS(file = "./HuVascAD_neuron.rds")
HuVascAD <- readRDS(file = "HuVascAD_astrocyte.rds")


pc_data <- data.frame(HuVascAD$treat, HuVascAD$gender,HuVascAD$sample, HuVascAD$age,HuVascAD$batch, HuVascAD$region, A[, "CDR"],A_varimax[, 1], A_varimax[, 2])
remove(HuVascAD)
colnames(pc_data) <- c("ident", "gender", "sample", "age", "batch","region", "CDR", "PC_1", "PC_2")
pc_data$gender <- as.factor(pc_data$gender)
pc_data$batch <- as.factor(pc_data$batch)
pc_data$sample <- as.factor(pc_data$sample)
pc_data$region <- as.factor(pc_data$region)
pc_data$ident <- as.factor(pc_data$ident)
pc_data$age <- as.numeric(pc_data$age)
pc_data$CDR <- as.numeric(pc_data$CDR)

subset <- which(pc_data$batch %in% c("1", "2", "3","4")) ### takes all samples
#subset <- which( pc_data$sample %in% c( "C1", "C7", "C8", "C9", "AD4", "AD5", "AD6", "AD7"))# pc_data$batch==4 & sample(n, size=1700) "AD4", "C7", "C8"

adjusted_counts <- B[, c(1,2,3,7)] %*% t(A[subset, c(1,2,3,7)]) + B_varimax[,1:2 ] %*% t(A_varimax[subset, 1:2])#B[, k] %*% t(A[subset, k]) #+ 
#adjusted_counts <- B_varimax[, -1] %*% t(A_varimax[subset,-1])
rownames(adjusted_counts) <- rownames(B)
colnames(adjusted_counts) <- rownames(A)[subset]

test_data_B <- CreateSeuratObject(counts = adjusted_counts )
remove(adjusted_counts)
test_data_B <- FindVariableFeatures(test_data_B)
test_data_B <- ScaleData(test_data_B)
test_data_B <- RunPCA(test_data_B)
test_data_B <- RunUMAP(test_data_B, dims = 1:30)
test_data_B$region <- pc_data$region[subset]
test_data_B$gender <- pc_data$gender[subset]
test_data_B$sample <- pc_data$sample[subset]
test_data_B$batch <- pc_data$batch[subset]
test_data_B$age <- pc_data$age[subset]
test_data_B$orig.ident <- pc_data$ident[subset]
test_data_B$age <- pc_data$age[subset]
test_data_B$CDR <- pc_data$CDR[subset]

DimPlot(test_data_B, dims = c(1,2), reduction = "umap", group.by= "sample") #, split.by = "region"


###################################################################
###################################################################
###################################################################
### Check the proportion of variance that can be explained by sample or region ###
rm(list=ls())
suppressPackageStartupMessages({
  library(Seurat)
  #library(knockoff)
  #library(MAST)
  library(Matrix)
  library(softImpute)
  library(tidyverse)
  library(RColorBrewer)
})

### Load the imputation/factor analysis results
### This is obtained by running the script "HuVascAD_downNULL_geneNULL_scsoftImpute_PC43_ll01_hippo_onlyFALSE_astrocyte_factors.R"
load("./factors/HuVascAD_downNULL_geneNULL_scsoftImpute_PC43_hippo_onlyFALSE_astrocyte_factors.rda")
A <- FinalResult[[1]]
B <- FinalResult[[2]]
n = nrow(A)
p = nrow(B)
r = ncol(B)
k = 7 # number of covariates used
PC = 43 # number of factors for astrocyte data

A_var <- apply(A, 2, var)
A_var
B_var <- apply(B, 2, var)
B_var

### Remove the covariates, extract the remaining latent factors
A_sub <- A[, (k+1):r]
B_sub <- B[, (k+1):r]

### Apply the varimax rotation
B_rotated <- varimax(B_sub)
B_varimax <- B_rotated$loadings
A_varimax <- A_sub %*% B_rotated$rotmat

colnames(B_varimax) <- paste0("Factor ", c(1:(r-k)))


datafilename = paste0("/Users/linxiliu/Dropbox/Research/Knockoff_scRNAseq/Data/seperate_data/", "HuVascAD_astrocyte.rds")
HuVascAD <- readRDS(file = datafilename)

meta_data <- data.frame(HuVascAD$treat, HuVascAD$gender,HuVascAD$sample, HuVascAD$age,HuVascAD$batch, HuVascAD$region, A[, "CDR"])
remove(HuVascAD)
colnames(meta_data) <- c("ident", "gender", "sample", "age", "batch","region", "CDR")
meta_data$gender <- as.numeric(meta_data$gender)
meta_data$batch <- as.factor(meta_data$batch)
meta_data$sample <- as.factor(meta_data$sample)
meta_data$region <- as.factor(meta_data$region)
meta_data$ident <- as.factor(meta_data$ident)
meta_data$age <- as.numeric(meta_data$age)
meta_data$CDR <- as.numeric(meta_data$CDR)

samp = model.matrix(y~x-1,data.frame(x = meta_data$sample, y = 1))
samp = samp[,-1]

batch = model.matrix(y~x-1,data.frame(x = meta_data$batch, y = 1))
batch = batch[,-1]

reg = model.matrix(y~x-1,data.frame(x = meta_data$region, y = 1))
reg = reg[,-1]

X_0 <- A[, 1:k]

### Fixed effect

library(lme4)

PC = 43
summary_out_var<- rep(0, PC)
summary_samp_var <- rep(0, PC)
summary_reg_var <- rep(0, PC)
summary_samp_var_random <- rep(0, PC)
summary_reg_var_random <- rep(0, PC)
for(k in 1:PC) {
  temp_data <- A_varimax[,1:k] %*% t(B_varimax[,1:k])
  B_0 <- solve(t(X_0)%*% X_0)%*% t(X_0)%*%temp_data
  temp_data_adj <- temp_data - X_0 %*% B_0
  temp_data_adj <- scale(temp_data_adj, center = TRUE, scale = FALSE)
  temp_data_svd <- svd.als(temp_data_adj, rank.max = k+2, lambda = 0, maxit = 200, thresh = 1e-08)
  out_var <- rep(0, length(temp_data_svd$d))
  samp_var <- rep(0, length(temp_data_svd$d))
  reg_var <- rep(0, length(temp_data_svd$d))
  samp_var_random <- rep(0, length(temp_data_svd$d))
  reg_var_random <- rep(0, length(temp_data_svd$d))
  for(l in 1:length(temp_data_svd$d)) {
    temp_outcome <- temp_data_svd$u[,l] * temp_data_svd$d[l]
    ### Fixed effects
    temp_frame_samp_fix <- data.frame(temp_outcome, samp)
    temp_frame_region_fix <- data.frame(temp_outcome, reg)
    temp_lm_samp <- lm(temp_outcome~., data=temp_frame_samp_fix)
    temp_lm_reg <- lm(temp_outcome~., data= temp_frame_region_fix)
    temp_total_var <- sum((temp_outcome - mean(temp_outcome))^2)
    temp_exp_var_samp <- temp_total_var*summary(temp_lm_samp)$r.squared
    temp_exp_var_reg <- temp_total_var*summary(temp_lm_reg)$r.squared
    samp_var[l] <- temp_exp_var_samp
    reg_var[l] <- temp_exp_var_reg
    out_var[l] <- temp_total_var
    
    ### random effects
    # sample
    temp_frame_samp_random <- data.frame(temp_outcome, meta_data$sample)
    temp_lmer_samp <- lmer(temp_outcome ~  1+(1 | meta_data.sample), data = temp_frame_samp_random)
    vc <- VarCorr(temp_lmer_samp)
    var_sample <- as.numeric(vc$meta_data.sample)
    var_resid <- sigma(temp_lmer_samp)^2
    samp_var_random[l] <- var_sample/(var_sample+var_resid)
    # region
    temp_frame_region_random <- data.frame(temp_outcome, meta_data$region)
    temp_lmer_reg <- lmer(temp_outcome ~ 1+  (1 | meta_data.region), data = temp_frame_region_random)
    vc <- VarCorr(temp_lmer_reg)
    var_reg <- as.numeric(vc$meta_data.reg)
    var_resid <- sigma(temp_lmer_reg)^2
    reg_var_random[l] <- var_reg/(var_reg+var_resid)
  }
  summary_out_var[k] <- sum(out_var)
  summary_samp_var[k] <- sum(samp_var)
  summary_reg_var[k] <- sum(reg_var)
  summary_samp_var_random[k] <- sum(out_var * samp_var_random)
  summary_reg_var_random[k] <- sum(out_var*reg_var_random)
}

FinalResult <- list(summary_out_var = summary_out_var,
                    summary_samp_var = summary_samp_var,
                    summary_reg_var = summary_reg_var,
                    summary_samp_var_random = summary_samp_var_random,
                    summary_reg_var_random = summary_reg_var_random)

outputname = paste0("Astrocyte_ANOVA.rda" )
save(FinalResult,file = outputname)

#load(file = "Astrocyte_ANOVA.rda" )