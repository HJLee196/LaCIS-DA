rm(list = ls())

# Set the working directory to the folder containing the .rda files
setwd("./")

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


list_rda = list.files(pattern = ".rda$", recursive = FALSE)

df_final_result = NULL

target_FDR_set = seq(from = 0.01, to = 0.3, length.out = 60)

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

list_rda = list.files(pattern = ".rda$", recursive = TRUE)
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



# cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#CC79A7", "#C0C0C0")
# q_val_label = c(
#   '0.005' = "q = 0.005",
#   '0.01' = "q = 0.01",
#   '0.05' = "q = 0.05",
#   '0.1' = "q = 0.1"
# )
# 
# q_val_test_label = c(
#   '0.005' = "q = 0.005",
#   '0.01' = "q = 0.01",
#   '0.05' = "q = 0.05",
#   '0.1' = "q = 0.1",
#   'MAST' = "MAST",
#   'wilcox' = "WRT",
#   'LR' = "LRT",
#   'LRT' = "LRT",
#   'LCD' = 'LCD',
#   'LCD1' = 'LCD',
#   'asdp' = "ASDP",
#   "decomp" = "Decomp",
#   "LR" = "LR",
#   "multidecomp" = "Multidecomp",
#   "multiLR" = "MultiLR",
#   "edecomp" = "Edecomp",
#   "eLR" = "ELR"
# )

df_final_result = df_final_result %>% select(-sgn)

df_final_result_all = rbind(df_final_result, df_final_result2)

df_final_result_all$FDP[is.na(df_final_result_all$FDP)] = 0

llam_temp = 
  df_final_result_all$llam %>%
  str_replace(pattern = "^llam", replacement = "") %>%
  as.numeric()

llam_temp[df_final_result_all$method == "BH"] = 0

df_final_result_all$llam = (llam_temp*25.2) %>% as.factor() #252.0678

# df_final_result$method = factor(df_final_result$method,
#                                 levels = c("asdp", "decomp", "LR","multidecomp", "multiLR", "edecomp", "eLR"))



# df_final_result$llam = (llam_temp*25.2) %>% as.factor() #252.0678

# df_final_result$FDP[is.na(df_final_result$FDP)] = 0

pd = position_dodge2(width = 0.7, preserve = "single", padding = 0)

font_size = 10

df_final_result_all %>% colnames()

summary_final_result_all = 
  df_final_result_all %>%
  filter(llam %in% c(0, 25.2, 75.6)) %>%
  group_by(method, test_stat, target_FDR, llam) %>%
  summarise(mean_FDP = mean(FDP), mean_Power = mean(Power)) %>%
  mutate(method_llam_test = str_c(method, llam, test_stat, sep = "_"))

summary_final_result_all$method %>% unique()
summary_final_result_all$llam %>% unique()

# scale_color_values_vector = 
#   c(
#     "BH_0_LR" = "red",
#     "decomp_25.2_LCD1" = "#0072B2", 
#     "decomp_50.4_LCD1" = "#004C7A",
#     "multidecomp_25.2_LCD1" = "#D4C200",
#     "multidecomp_50.4_LCD1" = "#BFA600",
#     "X = Y line" = "black")

# scale_color_values_vector <- 
#   c(
#     "BH_0_LR" = "red",
#     # blue: much larger lightness contrast
#     "decomp_25.2_LCD1" = "#56B4E9",
#     "decomp_50.4_LCD1" = "#003B73",
#     # green: lighter vs darker
#     "multidecomp_25.2_LCD1" = "#8DD3C7",
#     "multidecomp_50.4_LCD1" = "#1B7A6E",
#     "X = Y line" = "black"
#   ) 

scale_color_values_vector <- 
  c(
    "BH_0_LR" = "red",
    
    "decomp_25.2_LCD1" = "#56B4E9",
    "decomp_75.6_LCD1" = "#003B73",
    
    # clearer green, distinct from blue
    "multidecomp_25.2_LCD1" = "#7FBF7B",
    "multidecomp_75.6_LCD1" = "#238B45",
    
    "X = Y line" = "black"
  )

upower = 
  summary_final_result_all %>%
  filter(method != "edecomp") %>%
  ggplot(aes(x = target_FDR, y = mean_Power, colour = method_llam_test)) +
  geom_line() +
  scale_color_manual(
    values = scale_color_values_vector
    ) +
  theme(legend.position = "none",
        axis.text=element_text(size = font_size - 1),
        axis.title=element_text(size = font_size), #,face="bold"
        strip.text.x = element_text(size = font_size - 1)) +
  xlab(NULL) +
  ylab("Power") +
  xlim(c(0, 0.3))

bfdr =
  summary_final_result_all %>%
    filter(method != "edecomp") %>%
    ggplot(aes(x = target_FDR, y = mean_FDP, colour = method_llam_test)) +
    geom_line() +
    geom_line(
      data = data.frame(x = c(0, 0.3), y = c(0, 0.3)),
      aes(x = x, y = y, color = "X = Y line"),
      linetype = "dashed"
    ) +
    scale_color_manual(
      name = "Knockoff \nconstruction \n& test statistic",
      values = scale_color_values_vector,
      breaks = c(
        "BH_0_LR",
        "decomp_25.2_LCD1",
        "decomp_75.6_LCD1",
        "multidecomp_25.2_LCD1",
        "multidecomp_75.6_LCD1",
        "X = Y line"
      ),
      labels = c(
        "BH LRT",
        expression(decomp~(lambda == 25.2)~LCD),
        expression(decomp~(lambda == 75.6)~LCD),
        expression(multidecomp~(lambda == 25.2)~LCD),
        expression(multidecomp~(lambda == 75.6)~LCD),
        "Reference: X = Y"
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = font_size + 0.5),
      legend.title = element_text(size = font_size + 1.5),
      axis.text = element_text(size = font_size - 1),
      axis.title = element_text(size = font_size),
      strip.text.x = element_text(size = font_size - 1)
    ) +
    xlab("Target FDR") +
    ylab("Average FDP")


full_plot = plot_grid(upower,bfdr, labels = c("A","B"), align = "v", rel_heights = c(0.8,1),
                      nrow = 2,ncol = 1)

full_plot

jpeg(filename = "SimB2_Revision_260621_1.jpeg", width = 8, height = 6, units = "in",
     res = 300)
print(full_plot)
dev.off()


##########
