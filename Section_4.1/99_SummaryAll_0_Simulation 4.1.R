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

target_FDR_set = c(0.1, 0.01)

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

pd = position_dodge2(width = 0.7, preserve = "single", padding = 0)

font_size = 10

upower =
  df_final_result %>%
  filter((llam %in% c("25.2", "75.6"))) %>%
  ggplot(aes(x = llam, y = Power, fill = method)) +
  geom_boxplot(width = 0.7, position = pd) + #fatten = 0, 
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.7, linewidth = 0.2, position = pd, color = "#8B0000",
               aes(group = interaction(llam, method))) + #
  facet_grid(.~test_stat + target_FDR, labeller = as_labeller(q_val_test_label)) + 
  theme(legend.position = "none",
        axis.text=element_text(size = font_size - 1),
        axis.title=element_text(size = font_size), #,face="bold"
        strip.text.x = element_text(size = font_size - 1)) +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("") + 
  scale_fill_manual(breaks=c("asdp", "decomp", "LR","multidecomp", "multiLR", "edecomp", "eLR"),
                    values=cbPalette)

bfdr = 
  df_final_result %>%
  filter((llam %in% c("25.2", "75.6"))) %>%
  ggplot(aes(x = llam, y = FDP, fill = method)) +
  geom_boxplot(width = 0.7, position = pd) + #fatten = 0, 
  stat_summary(fun = mean, geom = "crossbar",
               width = 0.7, linewidth = 0.2, position = pd, color = "#8B0000",
               aes(group = interaction(llam, method))) + #
  facet_grid(.~test_stat + target_FDR, labeller = as_labeller(q_val_test_label)) + 
  theme(legend.position = "none",
        legend.title = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        axis.text=element_text(size=font_size-1),
        axis.title.y=element_text(size=font_size),
        axis.title.x=element_text(size=font_size+5), #,face="bold"
        strip.text.x = element_text(size = font_size - 1)) +
  coord_cartesian(ylim = c(0, 0.25)) +
  xlab(expression(lambda)) +
  theme(legend.position="bottom")+
  scale_fill_manual(breaks=c("asdp", "decomp", "LR","multidecomp", "multiLR", "edecomp", "eLR"),
                    values=cbPalette, name = "Knockoff construction")


full_plot = plot_grid(upower,bfdr, labels = c("A","B"), align = "v", rel_heights = c(1.4,1),
                      nrow = 2,ncol = 1)

full_plot

jpeg(filename = "Sim41_Figure1_251215.jpeg", width = 10, height = 8, units = "in",
     res = 300)
print(full_plot)
dev.off()

final_summarize = 
  df_final_result %>%
  group_by(test_stat, method, llam, target_FDR) %>%
  summarise(avg_FDP = mean(FDP),
            med_FDP = median(FDP),
            avg_Power = mean(Power))

write.csv(x = final_summarize, file = "Sim41_Figure1_251215.csv", row.names = F)
