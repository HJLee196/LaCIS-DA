rm(list = ls())

library(MASS)
library(tidyverse)
library(data.table)
library(cowplot)

loss <- function(proposal,l,alpha, n, p){
  
  L=matrix(data=0,nrow=10,ncol=p)
  for (i in 1:10){
    #wi=rWishart(n, p, diag(rgamma(p,proposal,1)))/n
    sd_candidate = matrix(sqrt(rgamma(p,proposal,proposal)), n, p, byrow = T)
    x1 = matrix(rnorm(n = n*p), n, p)*sd_candidate
    
    if (p<=n){
      l1=eigen(t(x1)%*%x1, symmetric = T, only.values = T)$values
      L[i,]=l1/n
    }
    
    if (p>n){
      l1=eigen(x1%*%t(x1), symmetric = T, only.values = T)$values
      L[i,1:n]=l1/n
    }  
    
  }
  l1=colMeans(L)
  
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  s1=lm(l[k]~l1[k]-1)$coef[[1]]
  l1=s1*l1
  return(sum((l1-l)[k]^2))
}

SQM_Gamma <- function(data, alpha = 0.2, beta_q = 0.1, n_mc = 500){
  
  data = as.matrix(data)
  
  n = nrow(data)
  p = ncol(data)
  
  S <- t(data) %*% data /n
  
  # system.time(
  #   l <- eigen(S)$values)  
  system.time(
    l <- eigen(S, symmetric = T, only.values = T)$values)
  
  # Computing Time
  start_time_opt = Sys.time()
  #
  o=optimize(loss,interval=c(0.1, 50),
             l=l,alpha=alpha,n=n,p=p)
  a=o$minimum
  
  # Computing Time
  end_time_opt = Sys.time()
  nsim_total_seconds = as.numeric(difftime(end_time_opt , start_time_opt, units = "secs"))
  nsim_minutes <- floor(nsim_total_seconds / 60)
  nsim_seconds <- round(nsim_total_seconds %% 60)
  print(sprintf("Spent Computation Time (Opt): %02d:%02d", nsim_minutes, nsim_seconds))
  #
  L=matrix(data=0,nrow=n_mc,ncol=p)
  
  # Computing Time
  start_time_nsim = Sys.time()
  #
  for (i in 1:n_mc){
    # x1=mvrnorm(n,rep(0,p),diag(rgamma(p,a,a)))
    sd_candidate = matrix(sqrt(rgamma(p,a,a)), n, p, byrow = T)
    x1 = matrix(rnorm(n = n*p), n, p)*sd_candidate
    
    if (p<=n){
      l1 = eigen(t(x1)%*%x1, symmetric = T, only.values = T)$values
      L[i,]=l1/n
    }
    if (p>n){
      l1 = eigen(x1%*%t(x1), symmetric = T, only.values = T)$values
      L[i,1:n]=l1/n
    } 
    
    if ((i %% 10) == 1){
      # Computing Time
      end_time_nsim = Sys.time()
      nsim_total_seconds = as.numeric(difftime(end_time_nsim , start_time_nsim, units = "secs"))
      nsim_minutes <- floor(nsim_total_seconds / 60)
      nsim_seconds <- round(nsim_total_seconds %% 60)
      cat(i, "th iteration", "\n")
      print(sprintf("Spent Computation Time: %02d:%02d", nsim_minutes, nsim_seconds))
      start_time_nsim = Sys.time()
      #
    }
  }
  
  l1=c()
  l2=c()
  l3=c()
  l4=c()
  for (i in 1:p){
    l1[i]=mean(L[,i])
    l2[i]=quantile(L[,i],1-beta_q)
    l3[i]=quantile(L[,i],beta_q)
    l4[i]=quantile(L[,i],0.5)
  }
  k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))
  s1=lm(l[k]~l1[k]-1)$coef[[1]]
  
  cat("Beta:", beta_q, "\n")
  cat("Estimated Rank:", sum(l>max(l2*s1)), "\n")
  
  L_max = apply(L, MARGIN = 1, FUN = max) 
  quantile_seq = seq(from = 0.05, to = 0.95, by = 0.05)
  Ls_quantile = quantile(L_max,probs = quantile_seq)*s1
  result_vec = 
    sapply(Ls_quantile, 
          FUN = function(Ls_temp, l_temp){
            sum(l_temp > Ls_temp)
          }, l_temp = l)
  result_mat = cbind(result_vec, names(result_vec))
  colnames(result_mat) = c("Rank", "Quantile")
  
  result_list = 
    list(L_mat = L,
         l_true = l,
         sigma2_est = s1,
         theta_est = a,
         result_mat = result_mat,
         alpha = alpha,
         beta_q = beta_q)
  
  return(result_list)
}

# Real Data
seed_num = 2017
set.seed(seed_num)

# Sim1
file_path = "./"
file_name = "Section52_xp_data_downNULL_regressed_new_p1.rda"
file_path_name = paste(file_path, file_name, sep = "/")

load(file = file_path_name)

real_mat = xp_data.matrix_regressed_new
rm(xp_data.matrix_regressed_new)

# Center the dataset
real_mat_nan = real_mat
real_mat_nan[real_mat == 0] = NaN
real_mat_mean =
  real_mat_nan %>%
  apply(2, FUN = mean, na.rm = T) %>%
  matrix(nrow = nrow(real_mat), ncol = ncol(real_mat), byrow = T)
real_mat_demean = real_mat - real_mat_mean
real_mat_demean[real_mat == 0] = 0
print("Check if the original matrix is centered well: ")
real_mat_mean %>% abs() %>% max() %>% print()
print("Check if the matrix is centered well: ")
real_mat_demean %>% apply(2, mean) %>% abs() %>% max() %>% print()
print("Data Dim:")
print(dim(real_mat_demean))

# Hyperparameter
alpha = 0.2
beta_q = 0.25
n_mc = 500

real_result_list = SQM_Gamma(data = real_mat_demean, alpha = alpha, beta_q = beta_q, n_mc = n_mc)

file_name = 
  paste("Result_list_alpha_", 
        str_pad(string = alpha*100, width = 3, side = "left", pad = "0"), 
        "_seed", seed_num,
        "_nmc", n_mc,
        sep = "")

write.csv(x = real_result_list$result_mat, file = paste(file_name, ".csv", sep = ""))
save(real_result_list, file = paste(file_name, ".rda", sep = ""))

cat("Estimate sigma2:", real_result_list$sigma2_est, "\n")
cat("Estimate theta:", real_result_list$theta_est, "\n")

########
# Plot #
########

## 1. QQplot

L_mat = real_result_list$L_mat
l_true = real_result_list$l_true
alpha = real_result_list$alpha
sigma2_est = real_result_list$sigma2_est
beta_q = real_result_list$beta_q
n = nrow(real_mat)
p = ncol(real_mat)

L_mean = L_mat %>% apply(MARGIN = 2, FUN = mean)
L_mean_scaled = L_mean*sigma2_est
k=floor(min(p,n)*alpha):floor(min(p,n)*(1-alpha))

df_qqplot = cbind(l_true[k], L_mean_scaled[k])
colnames(df_qqplot) = c("True", "Fitted")

result_qqplot =
  ggplot(df_qqplot, aes(x = True, y = Fitted)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(x = "Theoretical Eigenvalues",
       y = "Fitted Eigenvalues") +
  theme(legend.position = "none", 
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text=element_text(size=14),
        axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15), #,face="bold"
        strip.text.x = element_text(size = 14))

## 2. Threshold Plot
L_mat %>% dim()
Ls_mat = L_mat*sigma2_est
df_threshold = tibble(true_eigenvalues = l_true)
df_threshold$fitted_mean = Ls_mat %>% apply(2, mean)
df_threshold$fitted_high = Ls_mat %>% apply(2, quantile, probs = 1-beta_q)
df_threshold$fitted_low = Ls_mat %>% apply(2, quantile, probs = beta_q)
df_threshold$fitted_med = Ls_mat %>% apply(2, quantile, probs = 0.5)
df_threshold$ind = 1:nrow(df_threshold)

df_threshold = df_threshold %>% gather(key = "key", value = "value", -ind)

result_threshold_plot =
  df_threshold %>%
  subset(key != "true_eigenvalues") %>%
  ggplot(aes(x = ind, y = value, color = key)) +
  geom_line(linewidth = 0.5) +
  geom_point(data = subset(df_threshold, key == "true_eigenvalues"), 
             aes(color = key), size = 0.5) +
  scale_color_manual(
    values = c("fitted_mean" = "red", "fitted_high" = "orange", "fitted_low" = "blue", "fitted_med" = "green",
               "true_eigenvalues" = "black"),
    labels = c("Fitted Upper Quartile", "Fitted Lower Quartile", "Fitted Mean", "Fitted Median", "True Eigenvalues")   
  ) +
  geom_hline(yintercept = df_threshold %>% filter(key == "fitted_high") %>% dplyr::select(value) %>% max(), color = "orange", alpha = 0.5) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(),
        # legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        axis.text=element_text(size=14),
        axis.title.y=element_text(size=15),
        axis.title.x=element_text(size=15), #,face="bold"
        strip.text.x = element_text(size = 14)) +
  labs(y = "Eigenvalues",
       x = "Index") +
  guides(color = guide_legend(nrow = 2))

plot_name =
  paste("Result_plot",
        "_alpha",
        str_pad(string = alpha*100, width = 3, side = "left", pad = "0"), 
        "_seed", seed_num,
        "_nmc", n_mc,
        ".jpeg",
        sep = "")

full_plot = plot_grid(result_qqplot, result_threshold_plot, 
                      labels = c("A","B"), align = "v", rel_heights = c(1.4,1),
                      nrow = 2,ncol = 1)

full_plot

jpeg(filename = plot_name, width = 10, height = 8, units = "in",
     res = 300)
print(full_plot)
dev.off()
  


