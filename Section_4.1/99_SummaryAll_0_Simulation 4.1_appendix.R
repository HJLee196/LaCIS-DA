rm(list = ls())

# Set the working directory to the folder containing the .rda files
setwd("path/to/your/rda_folder")

library(tidyverse)
library(knockoff)
library(cowplot)

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


list_rda = list.files(pattern = ".rda$")

df_final_result = NULL

target_FDR_set = c(0.05)

# i_rda = 1
for (i_rda in 1:length(list_rda))
{
  rda_file_name_temp = list_rda[i_rda]
  
  method_name =
    rda_file_name_temp %>%
    str_extract(pattern = "method[[:alpha:]]*") %>%
    str_replace(pattern = "method", replacement = "")
  
  test_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "testUse[[:alpha:][:digit:]]*") %>%
    str_replace(pattern = "testUse", replacement = "")
  
  llam_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "llam[[:digit:].]*") %>%
    str_replace(pattern = "testUse", replacement = "")
  
  sgn_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "sgn[[:digit:]\\.]*") %>%
    str_replace(pattern = "sgn", replacement = "")
  
  if( method_name == "BH" | method_name == "naive")
  {
    next
  }
  
  load(rda_file_name_temp)
  
  # i_FinalResult = 1
  for (i_FinalResult in 1:length(FinalResult))
  {
    FinalResult_temp = FinalResult[[i_FinalResult]]
    
    W_imp_b_temp = FinalResult_temp$W_imp_b
    sig_coef_temp = FinalResult_temp$sig_coef
    
    # i_target_FDR = 1
    for (i_target_FDR in 1:length(target_FDR_set))
    {
      # col 
      # method # test stat # target_FDR # FDP # Power
      target_FDR_temp = target_FDR_set[i_target_FDR]
      
      # In the case of single knockoffs
      if (is.null(dim(W_imp_b_temp)))
      {
        threshold_temp <- knockoff.threshold(W_imp_b_temp, fdr = target_FDR_temp)
        
        selected_b_temp <- which(W_imp_b_temp >= threshold_temp)
        
        print((length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)) # FDP
        FDP_temp = (length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)
        
        print(sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)) # Power
        Power_temp = sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)
        
        df_sim_temp = data.frame(method = method_name, test_stat = test_name, sgn = sgn_name,
                                 target_FDR = target_FDR_temp, FDP = FDP_temp, Power = Power_temp,
                                 llam = llam_name)
        
        df_final_result = rbind(df_final_result, df_sim_temp)
      } else {
        # In the case of multiple knockoffs
        
        # 1. Multi # (1. Multi 2. ebh) 
        m_kos = dim(W_imp_b_temp)[1] - 1 # -1 for the orignal dataset (1 original 5 knockoffses)
        
        W_max_ind = apply(W_imp_b_temp,2,which.max)
        W_max = sapply(1:length(W_max_ind), function(x) W_imp_b_temp[W_max_ind[x],x])
        W_med = sapply(1:length(W_max_ind), function(x) median(W_imp_b_temp[-W_max_ind[x],x]))
        
        tau_stab <- (W_max - W_med)
        kap = (W_max_ind==1 & tau_stab>0)
        
        ts = sort(c(0, abs(tau_stab)))
        ratio = sapply(ts, function(t) (1/m_kos + 1/m_kos*sum((tau_stab >= t)*(1-kap)))/max(1, 
                                                                                            sum((tau_stab >= t)*(kap))))
        ok = which(ratio <= target_FDR_temp)
        threshold_temp <- ifelse(length(ok) > 0, ts[ok[1]], Inf)
        
        selected_b_temp <- which(tau_stab*kap >= threshold_temp)
        
        print((length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)) # FDR
        FDP_temp = (length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)
        
        print(sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)) # Power
        Power_temp = sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)
        
        df_sim_temp = data.frame(method = method_name, test_stat = test_name, sgn = sgn_name,
                                 target_FDR = target_FDR_temp, FDP = FDP_temp, Power = Power_temp,
                                 llam = llam_name)
        
        df_final_result = rbind(df_final_result, df_sim_temp)
        
        # 2. ebh # (1. Multi 2. ebh) 
        
        E_temp = matrix(NA, m_kos, dim(W_imp_b_temp)[2])
        for(m_ko in 1:m_kos){
          
          E_temp[m_ko,] = W_imp_b_temp[1,] - W_imp_b_temp[(m_ko+1),]
          tau = knockoff.threshold(E_temp[m_ko,],
                                   fdr = target_FDR_temp/2)
          E_temp[m_ko,] = (E_temp[m_ko,] >= tau) / (1 + sum(E_temp[m_ko,] <= -tau))
        }
        
        E_temp <- (dim(W_imp_b_temp)[2]) * colMeans(E_temp)
        
        selected_b_temp = ebh(E_temp, alpha = target_FDR_temp)$rej
        
        print((length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)) # FDR
        FDP_temp = (length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)
        
        print(sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)) # Power
        Power_temp = sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)
        
        method_name_ebh = 
          method_name %>% 
          str_replace(pattern = "^multi", replacement = "e")
        
        df_sim_temp = data.frame(method = method_name_ebh, test_stat = test_name, sgn = sgn_name,
                                 target_FDR = target_FDR_temp, FDP = FDP_temp, Power = Power_temp,
                                 llam = llam_name)
        
        df_final_result = rbind(df_final_result, df_sim_temp)
      }
    }
  }
}

cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#C0C0C0")
q_val_label = c(
  '0.005' = "q = 0.005",
  '0.01' = "q = 0.01",
  '0.05' = "q = 0.05",
  '0.1' = "q = 0.1"
)

q_val_test_label = c(
  '0.005' = "q = 0.005",
  '0.01' = "q = 0.01",
  '0.05' = "q = 0.05",
  '0.1' = "q = 0.1",
  'MAST' = "MAST",
  'wilcox' = "WRT",
  'LR' = "LRT",
  'LRT' = "LRT",
  'LCD' = 'LCD',
  'LCD1' = 'LCD',
  'asdp' = "ASDP",
  "decomp" = "Decomp", 
  "LR" = "LR",
  "multidecomp" = "Multidecomp", 
  "multiLR" = "MultiLR", 
  "edecomp" = "Edecomp", 
  "eLR" = "ELR"
)

df_final_result$method = factor(df_final_result$method,
                                levels = c("asdp", "decomp", "LR","multidecomp", "multiLR", "edecomp", "eLR"))

llam_temp = 
  df_final_result$llam %>%
  str_replace(pattern = "^llam", replacement = "") %>%
  as.numeric()

df_final_result$llam = (llam_temp*25.2) %>% as.factor()

df_final_result$FDP[is.na(df_final_result$FDP)] = 0

final_summarize = 
  df_final_result %>%
  group_by(test_stat, method, llam, target_FDR) %>%
  summarise(avg_FDP = mean(FDP),
            sd_FDP = sd(FDP),
            avg_Power = mean(Power),
            sd_Power = sd(Power))

col_values = c("avg_FDP", "sd_FDP", "avg_Power", "sd_Power")

final_summarize =
  final_summarize %>%
  pivot_wider(names_from = "test_stat",
              values_from = col_values,
              names_glue = "{test_stat}_{.value}")

col_idx = 
  (matrix(1:16, nrow = 4, ncol = 4, byrow = T) %>% 
     as.vector()) + 3

final_summarize_ordered = final_summarize[,c(1,2,3, col_idx)]

# write.csv(x = final_summarize_ordered, file = "Sim41_TableS4_asdp_smallLlam_251219.csv", row.names = F)

# Decomp Naive BH

list_rda = list.files(pattern = ".rda$")
idx_decompNaiveBH =
  list_rda %>%
  str_detect(pattern = "BH|naive")
list_rda = list_rda[idx_decompNaiveBH]

print(target_FDR_set)

df_final_result2 = NULL

for (i_rda in 1:length(list_rda))
{
  rda_file_name_temp = list_rda[i_rda]
  
  method_name =
    rda_file_name_temp %>%
    str_extract(pattern = "method[[:alpha:]]*") %>%
    str_replace(pattern = "method", replacement = "")
  
  test_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "testUse[[:alpha:][:digit:]]*") %>%
    str_replace(pattern = "testUse", replacement = "")
  
  llam_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "llam[[:digit:]]*") %>%
    str_replace(pattern = "testUse", replacement = "")
  
  load(rda_file_name_temp)
  
  
  # i_FinalResult = 1
  for (i_FinalResult in 1:length(FinalResult))
  {
    FinalResult_temp = FinalResult[[i_FinalResult]]
    
    if (method_name == "BH") ################# IMPORTANT ################ 
    {
      W_imp_b_temp = FinalResult_temp$W_imp
    } else {
      W_imp_b_temp = FinalResult_temp$W_imp_b
    }
    sig_coef_temp = FinalResult_temp$sig_coef
    
    # i_target_FDR = 1
    for (i_target_FDR in 1:length(target_FDR_set))
    {
      # col 
      # method # test stat # target_FDR # FDP # Power
      target_FDR_temp = target_FDR_set[i_target_FDR]
      
      if (method_name == "BH")
      {
        BH_p = p.adjust(p = W_imp_b_temp, method = "BH")
        selected_b_temp = (1:length(BH_p))[BH_p < target_FDR_temp]
        
        print((length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)) # FDR
        FDP_temp = (length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)
        
        print(sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)) # Power
        Power_temp = sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)
      } else {
        threshold_temp <- knockoff.threshold(W_imp_b_temp, fdr = target_FDR_temp)
        
        selected_b_temp <- which(W_imp_b_temp >= threshold_temp)
        
        print((length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)) # FDP
        FDP_temp = (length(selected_b_temp) - sum(selected_b_temp %in% sig_coef_temp))/length(selected_b_temp)
        
        print(sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)) # Power
        Power_temp = sum(selected_b_temp %in% sig_coef_temp)/length(sig_coef_temp)
      }
      
      df_sim_temp = data.frame(method = method_name, test_stat = test_name,
                               target_FDR = target_FDR_temp, FDP = FDP_temp, Power = Power_temp,
                               llam = llam_name)
      
      df_final_result2 = rbind(df_final_result2, df_sim_temp)
      
    }
  }
}


cbPalette2 <- c("#E69F00", "#F0E442", "#009E73")

q_val_test_label = c(
  '0.005' = "q = 0.005",
  '0.01' = "q = 0.01",
  '0.05' = "q = 0.05",
  '0.1' = "q = 0.1",
  'MAST' = "MAST",
  'wilcox' = "WRT",
  'LR' = "LRT",
  'LRT' = "LRT",
  'LCD' = 'NOLCD',
  "LCD1" = "LCD",
  '1' = "sgn = 1",
  '1.5' = "sgn = 1.5",
  '2' = "sgn = 2",
  '2.5' = "sgn = 2.5",
  '3' = "sgn = 3"
)

df_final_result2$method = 
  df_final_result2$method %>%
  str_replace_all(pattern = "naive", replacement = "Gaussian")

df_final_result2$method = factor(df_final_result2$method,
                                levels = c("Gaussian", "BH"))

llam_temp = 
  df_final_result2$llam %>%
  str_replace(pattern = "^llam", replacement = "") %>%
  as.numeric()

df_final_result2$llam = (llam_temp*25.2) %>% as.factor()

df_final_result2$FDP[is.na(df_final_result2$FDP)] = 0

# How well LCD is doing
final_summarize2 = 
  df_final_result2 %>%
  # filter(test_stat == "LCD") %>%
  group_by(test_stat, method, llam, target_FDR) %>%
  summarise(avg_FDP = mean(FDP),
            sd_FDP = sd(FDP),
            avg_Power = mean(Power),
            sd_Power = sd(Power))

final_summarize2 =
  final_summarize2 %>%
  pivot_wider(names_from = "test_stat",
              values_from = col_values,
              names_glue = "{test_stat}_{.value}")

col_idx = 
  (matrix(1:16, nrow = 4, ncol = 4, byrow = T) %>% 
     as.vector()) + 3

final_summarize2_ordered = final_summarize2[,c(1,2,3, col_idx)]
#Check if they have the same values for different llams

final_summarize2_ordered = 
  final_summarize2_ordered %>%
  filter(llam == unique(final_summarize2_ordered$llam)[1]) %>%
  mutate(llam = as.factor(0))

final_summarize_ordered_full =
  rbind(final_summarize_ordered,
        final_summarize2_ordered)

write.csv(x = final_summarize_ordered_full, 
          file = "Sim41_TableS4_251215.csv", row.names = F)


