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


list_rda = list.files(pattern = ".rda$", recursive = T)

df_final_result = NULL

# target_FDR_set = c(0.1, 0.01)

# i_rda = 1
for (i_rda in 1:length(list_rda))
{
  rda_file_path_temp = list_rda[i_rda]
  
  rda_file_name_temp = 
    rda_file_path_temp %>% 
    str_extract(pattern = "[^///]*.rda$")
  
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
  
  cell_name =
    rda_file_name_temp %>%
    str_extract(pattern = "cell[[:alpha:]]*") %>%
    str_replace(pattern = "cell", replacement = "")
  
  PC_name =
    rda_file_name_temp %>%
    str_extract(pattern = "(PC)([[:digit:]]*)", group = 2)
  
  cell_PC_name = 
    str_c(cell_name, PC_name, sep = "_")
  
  if( method_name == "BH" | method_name == "naive")
  {
    next
  }
  
  load(rda_file_path_temp)
  
  selected = FinalResult[[1]]$selected
  selected_b = FinalResult[[1]]$selected_b
  selected_b_ebh = FinalResult[[1]]$selected_b_ebh
  
  df_sim_temp = data.frame(method = method_name, 
                           test_stat = test_name, cell_PC_name = cell_PC_name,
                           target_FDR = 0.1, num_DEG = length(selected))
  
  df_final_result = rbind(df_final_result, df_sim_temp)
  
  df_sim_temp = data.frame(method = paste0(method_name, "_b"), 
                           test_stat = test_name, cell_PC_name = cell_PC_name,
                           target_FDR = 0.1, num_DEG = length(selected_b))
  
  df_final_result = rbind(df_final_result, df_sim_temp)
  
  df_sim_temp = data.frame(method = method_name %>% str_replace(pattern = "^multi", replacement = "e-"), 
                           test_stat = test_name, cell_PC_name = cell_PC_name,
                           target_FDR = 0.1, num_DEG = length(selected_b_ebh))
  
  df_final_result = rbind(df_final_result, df_sim_temp)
  
}

cell_PC_level_set = NULL
cell_level_set = c("astrocyte", "microglia", "neuron", "oligo", "pericyte")
unique_cell_PC_name = df_final_result$cell_PC_name %>% unique()
for (i_cell in 1:length(cell_level_set)){
  cell_temp = cell_level_set[i_cell]
  idx_cell_PC_temp = 
    unique_cell_PC_name %>% 
    str_detect(pattern = str_c("^", cell_temp))
  cell_PC_temp = unique_cell_PC_name[idx_cell_PC_temp]
  cell_PC_level_set = c(cell_PC_level_set, cell_PC_temp)
}


df_final_result$cell_PC_name =
  factor(df_final_result$cell_PC_name,
         levels = cell_PC_level_set)

df_final_result$test_stat = 
  factor(df_final_result$test_stat,
         levels = c("LCD1", "LR", "MAST", "wilcox"))

df_final_result$method = 
  factor(df_final_result$method,
         levels = c("multidecomp", "multidecomp_b", "multiLR", "multiLR_b", "e-decomp", "e-LR"))

final_summarize = 
  df_final_result %>%
  select(-target_FDR) %>%
  filter(method != "multiLR") %>%
  filter(method != "multidecomp") %>%
  spread(key = cell_PC_name, value = num_DEG)

final_summarize

write.csv(x = final_summarize, file = "Sim52_Table1.csv", row.names = F)
