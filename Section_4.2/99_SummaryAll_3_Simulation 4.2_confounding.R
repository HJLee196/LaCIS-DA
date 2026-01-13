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

target_FDR_set = c(0.1)

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
    str_extract(pattern = "llam[[:digit:]]*") %>%
    str_replace(pattern = "testUse", replacement = "")
  
  sgn_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "sgn[[:digit:]\\.]*") %>%
    str_replace(pattern = "sgn", replacement = "")
  
  down_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "down[[:digit:][:alpha:]\\.]*") %>%
    str_replace(pattern = "down", replacement = "")
  
  batch_name =
    rda_file_name_temp %>%
    str_extract(pattern = "batch[[:alpha:]\\.]*") %>%
    str_replace(pattern = "batch", replacement = "")
  
  impCDR_name =
    rda_file_name_temp %>%
    str_extract(pattern = "impCDR[[:alpha:]\\.]*") %>%
    str_replace(pattern = "impCDR", replacement = "")
  
  compCDR_name =
    rda_file_name_temp %>%
    str_extract(pattern = "compCDR[[:alpha:]\\.]*") %>%
    str_replace(pattern = "compCDR", replacement = "")

  
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
                                 llam = llam_name, down = down_name, 
                                 batch = batch_name, impCDR = impCDR_name, compCDR = compCDR_name)
        
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
                                 llam = llam_name, down = down_name, 
                                 batch = batch_name, impCDR = impCDR_name, compCDR = compCDR_name)
        
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
                                 llam = llam_name, down = down_name, 
                                 batch = batch_name, impCDR = impCDR_name, compCDR = compCDR_name)
        
        df_final_result = rbind(df_final_result, df_sim_temp)
      }
    }
  }
}

df_final_result$method = factor(df_final_result$method,
                                levels = c("asdp", "decomp", "LR","multidecomp", "multiLR", "edecomp", "eLR"))
df_final_result$llam = factor(df_final_result$llam,
                              levels = c("llam0", "llam1", "llam2", "llam3"))
df_final_result$FDP[is.na(df_final_result$FDP)] = 0

final_summarize_type1 = 
  df_final_result %>%
  group_by(llam, batch, impCDR, compCDR) %>%
  summarise(avg_FDP = mean(FDP),
            med_FDP = median(FDP),
            sd_FDP = sd(FDP),
            avg_Power = mean(Power),
            sd_Power = sd(Power))


final_summarize_name = 
  str_c("Section42_Table1_LR_nSim20",
        "_Down", down_name)

###########
# AVG FDP #
###########

final_summarize_type2_FDP =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(batch_llam = str_c(batch, llam, sep = "_"),
         imp_comp = str_c(impCDR, compCDR, sep = "_")) %>%
  mutate(imp_comp = imp_comp %>%
           str_replace(pattern = "^without_without",
                       replacement = "1_without_without") %>%
           str_replace(pattern = "^without_with",
                       replacement = "2_without_with") %>%
           str_replace(pattern = "^with_with",
                       replacement = "3_with_with")
           ) %>%
  dplyr::select(avg_FDP, batch_llam, imp_comp) %>%
  # mutate(avg_FDP = round(avg_FDP, 2)) %>%
  spread(key = imp_comp, value = avg_FDP) %>%
  mutate(batch_llam_sort = batch_llam %>%
           str_replace(pattern = "^permute_",
                       replacement = "2_permute_") %>%
           str_replace(pattern = "^with_",
                       replacement = "1_with_") %>%
           str_replace(pattern = "^without_",
                       replacement = "3_without")) %>%
  arrange(batch_llam_sort) %>%
  dplyr::select(-batch_llam_sort)

write.csv(x = final_summarize_type2_FDP, 
          file = str_c(final_summarize_name,
                       "_Avg_FDP", 
                       ".csv"), row.names = F)
  

###########
# SD FDP #
###########

final_summarize_type2_FDPSD =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(batch_llam = str_c(batch, llam, sep = "_"),
         imp_comp = str_c(impCDR, compCDR, sep = "_")) %>%
  mutate(imp_comp = imp_comp %>%
           str_replace(pattern = "^without_without",
                       replacement = "1_without_without") %>%
           str_replace(pattern = "^without_with",
                       replacement = "2_without_with") %>%
           str_replace(pattern = "^with_with",
                       replacement = "3_with_with")
  ) %>%
  dplyr::select(sd_FDP, batch_llam, imp_comp) %>%
  # mutate(sd_FDP = round(sd_FDP, 2)) %>%
  spread(key = imp_comp, value = sd_FDP) %>%
  mutate(batch_llam_sort = batch_llam %>%
           str_replace(pattern = "^permute_",
                       replacement = "2_permute_") %>%
           str_replace(pattern = "^with_",
                       replacement = "1_with_") %>%
           str_replace(pattern = "^without_",
                       replacement = "3_without")) %>%
  arrange(batch_llam_sort) %>%
  dplyr::select(-batch_llam_sort)

write.csv(x = final_summarize_type2_FDPSD, 
          file = str_c(final_summarize_name,
                       "_SD_FDP", 
                       ".csv"), row.names = F)

###########
# AVG Power #
###########

final_summarize_type2_Power =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(batch_llam = str_c(batch, llam, sep = "_"),
         imp_comp = str_c(impCDR, compCDR, sep = "_")) %>%
  mutate(imp_comp = imp_comp %>%
           str_replace(pattern = "^without_without",
                       replacement = "1_without_without") %>%
           str_replace(pattern = "^without_with",
                       replacement = "2_without_with") %>%
           str_replace(pattern = "^with_with",
                       replacement = "3_with_with")
  ) %>%
  dplyr::select(avg_Power, batch_llam, imp_comp) %>%
  # mutate(avg_Power = round(avg_Power, 2)) %>%
  spread(key = imp_comp, value = avg_Power) %>%
  mutate(batch_llam_sort = batch_llam %>%
           str_replace(pattern = "^permute_",
                       replacement = "2_permute_") %>%
           str_replace(pattern = "^with_",
                       replacement = "1_with_") %>%
           str_replace(pattern = "^without_",
                       replacement = "3_without")) %>%
  arrange(batch_llam_sort) %>%
  dplyr::select(-batch_llam_sort)

write.csv(x = final_summarize_type2_Power, 
          file = str_c(final_summarize_name,
                       "_Avg_Power", 
                       ".csv"), row.names = F)

###########
# SD Power #
###########

final_summarize_type2_PowerSD =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(batch_llam = str_c(batch, llam, sep = "_"),
         imp_comp = str_c(impCDR, compCDR, sep = "_")) %>%
  mutate(imp_comp = imp_comp %>%
           str_replace(pattern = "^without_without",
                       replacement = "1_without_without") %>%
           str_replace(pattern = "^without_with",
                       replacement = "2_without_with") %>%
           str_replace(pattern = "^with_with",
                       replacement = "3_with_with")
  ) %>%
  dplyr::select(sd_Power, batch_llam, imp_comp) %>%
  # mutate(sd_Power = round(sd_Power, 2)) %>%
  spread(key = imp_comp, value = sd_Power) %>%
  mutate(batch_llam_sort = batch_llam %>%
           str_replace(pattern = "^permute_",
                       replacement = "2_permute_") %>%
           str_replace(pattern = "^with_",
                       replacement = "1_with_") %>%
           str_replace(pattern = "^without_",
                       replacement = "3_without")) %>%
  arrange(batch_llam_sort) %>%
  dplyr::select(-batch_llam_sort)

write.csv(x = final_summarize_type2_PowerSD, 
          file = str_c(final_summarize_name,
                       "_SD_Power", 
                       ".csv"), row.names = F)



#############
# BH 
#############

rm(list = ls())

# Set the working directory to the folder containing the .rda files from the BH procedure
setwd("path/to/your/rda_folder")

list_rda = list.files(pattern = ".rda$")

df_final_result2 = NULL

target_FDR_set = c(0.1)

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
    str_extract(pattern = "llam[[:digit:]]*") %>%
    str_replace(pattern = "testUse", replacement = "")
  
  sgn_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "sgn[[:digit:]\\.]*") %>%
    str_replace(pattern = "sgn", replacement = "")
  
  down_name = 
    rda_file_name_temp %>%
    str_extract(pattern = "down[[:digit:][:alpha:]\\.]*") %>%
    str_replace(pattern = "down", replacement = "")
  
  batch_name =
    rda_file_name_temp %>%
    str_extract(pattern = "batch[[:alpha:]\\.]*") %>%
    str_replace(pattern = "batch", replacement = "")
  
  impCDR_name =
    rda_file_name_temp %>%
    str_extract(pattern = "impCDR[[:alpha:]\\.]*") %>%
    str_replace(pattern = "impCDR", replacement = "")
  
  compCDR_name =
    rda_file_name_temp %>%
    str_extract(pattern = "compCDR[[:alpha:]\\.]*") %>%
    str_replace(pattern = "compCDR", replacement = "")
  
  load(rda_file_name_temp)
  
  # i_FinalResult = 1
  for (i_FinalResult in 1:length(FinalResult))
  {
    FinalResult_temp = FinalResult[[i_FinalResult]]
    
    if (method_name == "BH")
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
                               llam = llam_name, down = down_name, 
                               batch = batch_name, impCDR = impCDR_name, compCDR = compCDR_name)
      
      df_final_result2 = rbind(df_final_result2, df_sim_temp)
      
    }
  }
}


df_final_result2$method = factor(df_final_result2$method,
                                levels = c("asdp", "decomp", "LR","multidecomp", "multiLR", "edecomp", "eLR", "BH"))
df_final_result2$llam = factor(df_final_result2$llam,
                               levels = c("llam0", "llam1", "llam2", "llam3"))
df_final_result2$FDP[is.na(df_final_result2$FDP)] = 0

final_summarize_type1 = 
  df_final_result2 %>%
  filter(llam == "llam0") %>%
  group_by(test_stat, batch, impCDR, compCDR) %>%
  summarise(avg_FDP = mean(FDP),
            med_FDP = median(FDP),
            sd_FDP = sd(FDP),
            avg_Power = mean(Power),
            sd_Power = sd(Power))

final_summarize_name = 
  str_c("Section42_Table2_BH",
        "_Down", down_name)

###########
# AVG FDP #
###########

final_summarize_type2_FDP =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(test_batch = str_c(test_stat, batch, sep = "_")) %>%
  mutate(compCDR = compCDR %>%
         str_replace(pattern = "^with$",
                     replacement = "2_with") %>%
         str_replace(pattern = "^without$",
                     replacement = "1_without")) %>%
  dplyr::select(avg_FDP, test_batch, compCDR) %>%
  # mutate(avg_FDP = round(avg_FDP, 4)) %>%
  spread(key = compCDR, value = avg_FDP) %>%
  mutate(test_batch_sort = test_batch %>%
         str_replace(pattern = "_permute$",
                     replacement = "_2_permute") %>%
         str_replace(pattern = "_with$",
                     replacement = "_1_with") %>%
         str_replace(pattern = "_without$",
                     replacement = "_3_without")) %>%
  arrange(test_batch_sort) %>%
  dplyr::select(-test_batch_sort)

write.csv(x = final_summarize_type2_FDP, 
          file = str_c(final_summarize_name,
                       "_Avg_FDP", 
                       ".csv"), row.names = F)

###########
# SD FDP #
###########

final_summarize_type2_FDPSD =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(test_batch = str_c(test_stat, batch, sep = "_")) %>%
  mutate(compCDR = compCDR %>%
           str_replace(pattern = "^with$",
                       replacement = "2_with") %>%
           str_replace(pattern = "^without$",
                       replacement = "1_without")) %>%
  dplyr::select(sd_FDP, test_batch, compCDR) %>%
  # mutate(sd_FDP = round(sd_FDP, 4)) %>%
  spread(key = compCDR, value = sd_FDP) %>%
  mutate(test_batch_sort = test_batch %>%
           str_replace(pattern = "_permute$",
                       replacement = "_2_permute") %>%
           str_replace(pattern = "_with$",
                       replacement = "_1_with") %>%
           str_replace(pattern = "_without$",
                       replacement = "_3_without")) %>%
  arrange(test_batch_sort) %>%
  dplyr::select(-test_batch_sort)

write.csv(x = final_summarize_type2_FDPSD, 
          file = str_c(final_summarize_name,
                       "_SD_FDP", 
                       ".csv"), row.names = F)

###########
# AVG Power #
###########

final_summarize_type2_Power =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(test_batch = str_c(test_stat, batch, sep = "_")) %>%
  mutate(compCDR = compCDR %>%
           str_replace(pattern = "^with$",
                       replacement = "2_with") %>%
           str_replace(pattern = "^without$",
                       replacement = "1_without")) %>%
  dplyr::select(avg_Power, test_batch, compCDR) %>%
  # mutate(avg_Power = round(avg_Power, 4)) %>%
  spread(key = compCDR, value = avg_Power) %>%
  mutate(test_batch_sort = test_batch %>%
           str_replace(pattern = "_permute$",
                       replacement = "_2_permute") %>%
           str_replace(pattern = "_with$",
                       replacement = "_1_with") %>%
           str_replace(pattern = "_without$",
                       replacement = "_3_without")) %>%
  arrange(test_batch_sort) %>%
  dplyr::select(-test_batch_sort)

write.csv(x = final_summarize_type2_Power, 
          file = str_c(final_summarize_name,
                       "_Avg_Power", 
                       ".csv"), row.names = F)

###########
# SD Power #
###########

final_summarize_type2_PowerSD =
  final_summarize_type1 %>%
  ungroup() %>%
  mutate(test_batch = str_c(test_stat, batch, sep = "_")) %>%
  mutate(compCDR = compCDR %>%
           str_replace(pattern = "^with$",
                       replacement = "2_with") %>%
           str_replace(pattern = "^without$",
                       replacement = "1_without")) %>%
  dplyr::select(sd_Power, test_batch, compCDR) %>%
  # mutate(sd_Power = round(sd_Power, 4)) %>%
  spread(key = compCDR, value = sd_Power) %>%
  mutate(test_batch_sort = test_batch %>%
           str_replace(pattern = "_permute$",
                       replacement = "_2_permute") %>%
           str_replace(pattern = "_with$",
                       replacement = "_1_with") %>%
           str_replace(pattern = "_without$",
                       replacement = "_3_without")) %>%
  arrange(test_batch_sort) %>%
  dplyr::select(-test_batch_sort)

write.csv(x = final_summarize_type2_PowerSD, 
          file = str_c(final_summarize_name,
                       "_SD_Power", 
                       ".csv"), row.names = F)