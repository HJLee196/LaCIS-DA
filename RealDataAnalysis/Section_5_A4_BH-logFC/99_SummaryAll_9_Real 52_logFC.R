rm(list = ls())

library(tidyverse)
library(reshape2)
library(knockoff)

# Set the working directory to the folder containing the .rda files
setwd("path/to/your/rda_folder")
source("UtilityFunctions.R")

lfc = c(0, 0.05, 0.1, 0.25, 0.5)
q_set = c(0.1, 0.05)

list_rda = list.files(pattern = ".rda$", recursive = T)

selection_main_all = NULL

for(i_rda in 1:length(list_rda))
{
  rda_file_path_temp = list_rda[i_rda]
  
  load(file = rda_file_path_temp)
  
  rda_file_name_temp = 
    rda_file_path_temp %>% 
    str_extract(pattern = "[^///]*.rda$")
  
  cell_name =
    rda_file_name_temp %>%
    str_extract(pattern = "cell[[:alpha:]]*") %>%
    str_replace(pattern = "cell", replacement = "")
  
  PC_name =
    rda_file_name_temp %>%
    str_extract(pattern = "(PC)([[:digit:]]*)", group = 2)
  
  selection_main_temp = 
    data.frame(cell_name = paste(cell_name, PC_name, sep = "_"),
               lfc = lfc)
  
  selection_mast = matrix(0, nrow = length(lfc), ncol = length(q_set))
  selection_lr = matrix(0, nrow = length(lfc), ncol = length(q_set))
  selection_wrt = matrix(0, nrow = length(lfc), ncol = length(q_set))
  num_genes_mast = matrix(0, nrow = length(lfc), ncol = length(q_set))
  num_genes_lr = matrix(0, nrow = length(lfc), ncol = length(q_set))
  num_genes_wrt = matrix(0, nrow = length(lfc), ncol = length(q_set))
  
  colnames(selection_lr) = rep("lr", length(q_set))
  colnames(selection_mast) = rep("mast", length(q_set))
  colnames(selection_wrt) = rep("wrt", length(q_set))
  
  for(i_q in 1:length(q_set))
  {
    for(i_lfc in 1:length(lfc)) 
    {
      q_temp_value = q_set[i_q]
      lfc_value = lfc[i_lfc]
      
      # MAST
      q_value_matrix = FinalResult[[1]]$q_value_mast
      all_lfc = as.numeric(rownames(q_value_matrix))
      row_ind <- which(all_lfc == lfc_value)
      selection_ind <- which(q_value_matrix[row_ind, ] < q_temp_value)
      
      selection_mast[i_lfc, i_q] <- length(selection_ind)
      num_genes_mast[i_lfc, i_q] = sum(q_value_matrix[row_ind, ] <2)
      
      # LR
      q_value_matrix = FinalResult[[1]]$q_value_lr
      all_lfc = as.numeric(rownames(q_value_matrix))
      row_ind <- which(all_lfc == lfc_value)
      selection_ind <- which(q_value_matrix[row_ind, ] < q_temp_value)
      
      selection_lr[i_lfc, i_q] <- length(selection_ind)
      num_genes_lr[i_lfc, i_q] = sum(q_value_matrix[row_ind, ] <2)
      
      # WRT
      q_value_matrix = FinalResult[[1]]$q_value_wrt
      all_lfc = as.numeric(rownames(q_value_matrix))
      row_ind <- which(all_lfc == lfc_value)
      selection_ind <- which(q_value_matrix[row_ind, ] < q_temp_value)
      
      selection_wrt[i_lfc, i_q] <- length(selection_ind)
      num_genes_wrt[i_lfc, i_q] = sum(q_value_matrix[row_ind, ] <2)
      
    }
    colnames(selection_lr)[i_q] = paste(colnames(selection_lr)[i_q], q_temp_value, sep = "_")
    colnames(selection_mast)[i_q] = paste(colnames(selection_mast)[i_q], q_temp_value, sep = "_")
    colnames(selection_wrt)[i_q] = paste(colnames(selection_wrt)[i_q], q_temp_value, sep = "_")
  }
  # Number of Genes
  check_num_genes_mast_lr = num_genes_mast - num_genes_lr
  check_num_genes_mast_wrt = num_genes_mast - num_genes_wrt
  sum_err = sum(abs(check_num_genes_mast_lr) + abs(check_num_genes_mast_wrt))
  
  if (sum_err != 0){
    print("################## ERROR ##################")
  } else {
    num_of_genes = num_genes_mast[,1, drop = F]
    colnames(num_of_genes) = "num_of_genes"
  }
  
  selection_all = 
    do.call(cbind, 
            list(num_of_genes,
                 selection_lr, 
                 selection_mast, 
                 selection_wrt))
  
  selection_main_temp = 
    cbind(selection_main_temp, selection_all)
  
  selection_main_all = rbind(selection_main_all, selection_main_temp)
}

selection_main_all

write.csv(x = selection_main_all, file = "Sim5_logFC_TableS3.csv", row.names = F)
