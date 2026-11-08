rm(list = ls())

library(tidyverse)

setwd("./")

# Read R script 48
mother_file_name = "48_mother_code_BH_realdata.R"

Rcodetext <- readLines(mother_file_name)

cell_name_set = c("astrocyte", "oligo", "pericyte", "microglia", "neuron")
PC_set = c("43", "64", "76", "11", "13")

llam_set = 0.1 %>% as.character() #c(0, 0.1, 0.2, 0.3)

nSim = 1
down = NULL
seed_num = 17
# file_path = "path/to/your/rds_folder"
file_path = "./"
max_iter_imp = 100

where_seed_num = grep("^seed_num[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_down = grep("^down[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_gene.index = grep("^gene.index[[:blank:]]*=[[:blank:]]*",Rcodetext)
# where_sign_strength = grep("^sign_strength[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_n_signal = grep("^n_signal[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_PC = grep("^PC[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_llam = grep("^llam[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_nSim = grep("^nSim[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_testUse = grep("^testUse[[:blank:]]*=[[:blank:]]*",Rcodetext)
where_max_iter_imp = grep("^max_iter_imp[[:blank:]]*=[[:blank:]]*",Rcodetext)

# where_libpath = grep("\\.libPaths",Rcodetext)
where_readRDS = 
  grep("^HuVascAD = readRDS", Rcodetext)

where_file_path = 
  grep("^file_path", Rcodetext)
where_cell_name =
  grep("^cell_name", Rcodetext)

# Change Single
Rcodetext[where_nSim] =
  Rcodetext[where_nSim] %>%
  str_replace(pattern = "^nSim[[:blank:]]*=[[:blank:]]*[[:digit:]]*", 
              replacement = paste0("nSim = ", nSim))

# Rcodetext[where_llam] =
#   Rcodetext[where_llam] %>%
#   str_replace(pattern = "^llam[[:blank:]]*=[[:blank:]]*[[:digit:]\\.]*", replacement = "llam = 0.1")

Rcodetext[where_down] =
  Rcodetext[where_down] %>%
  str_replace(pattern = "(^down[[:blank:]]*=[[:blank:]]*)([[:alnum:]]*)",
              replacement = paste0("\\1", ifelse(is.null(down), yes = "NULL", no = down)))

# Rcodetext[where_libpath] =
#   Rcodetext[where_libpath] %>%
#   str_replace("^#[[:blank:]]*", replacement = "")

Rcodetext[where_file_path] = 
  paste0("file_path = '",
         file_path,
         "'")

Rcodetext[where_seed_num] =
  Rcodetext[where_seed_num] %>%
  str_replace(pattern = "^seed_num[[:blank:]]*=[[:blank:]]*[[:digit:]]*", 
              replacement = paste0("seed_num = ", seed_num))

Rcodetext[where_max_iter_imp] = 
  Rcodetext[where_max_iter_imp] %>%
  str_replace(pattern = "(^max_iter_imp[[:blank:]]*=[[:blank:]]*)([[:digit:]]*)", 
              replacement = str_c("\\1", max_iter_imp))

# Read sh script
mother_file_name_sh =
  str_replace(mother_file_name, pattern = "\\.R", replacement = "\\.sh")
shcodetext <- readLines(mother_file_name_sh)

where_Rscript = grep("^Rscript",shcodetext)
where_jobName = grep("job-name",shcodetext)
where_mem = grep("mem", shcodetext)
where_time = grep("time", shcodetext)

# Change Single
shcodetext[where_mem] = 
  shcodetext[where_mem] %>%
  str_replace(pattern = "mem=[[:digit:]]*g",
              replacement = "mem=256g")

shcodetext[where_time] =
  shcodetext[where_time] %>%
  str_replace(pattern = "[[:digit:]]*:[[:digit:]]*:[[:digit:]]*",
              replacement = "23:00:00")


# Modify the parameters and save it.

for (i_llam in 1:length(llam_set))
{
for (i_PC_cell_name in 1:length(cell_name_set))
{
  # Initialize
  Rcodetext_child = Rcodetext
  shcodetext_child = shcodetext
  
  folder_name = 
    mother_file_name %>%
    str_to_lower() %>%
    str_replace(pattern = ".r$", replacement = "")
  
  if (!dir.exists(folder_name)){
    dir.create(folder_name)
  }
  
  # Change Parameters for the R script
  
  # llam
  llam_temp = llam_set[i_llam]
  Rcodetext_child[where_llam] =
    Rcodetext_child[where_llam] %>% 
    str_replace(pattern = "(=[[:blank:]]*)([[:digit:]\\.]*)", replacement = str_c("\\1", llam_temp))

  # PC and Cell_name
  cell_name_temp = cell_name_set[i_PC_cell_name]
  PC_temp = PC_set[i_PC_cell_name]
  Rcodetext_child[where_PC] = 
    Rcodetext_child[where_PC] %>%
    str_replace(pattern = "^(PC[[:blank:]]*=[[:blank:]]*)([[:digit:]\\.]*)",
                replacement = str_c("\\1", PC_temp))
  
  Rcodetext_child[where_cell_name] =
    Rcodetext[where_cell_name] %>%
    str_replace(pattern = "^(cell_name[[:blank:]]*=[[:blank:]]*)([[:alpha:][:punct:]]*)", 
                replacement = str_c("\\1", '\\\"', cell_name_temp, '\\\"'))

  
  # Save the R file.
  file_name_R =
    str_replace(mother_file_name, pattern = "mother", replacement = "child") %>%
    str_replace(pattern = "\\.R", replacement = str_c("_cell", cell_name_temp, "\\.R")) %>%
    str_replace(pattern = "\\.R", replacement = str_c("_down", ifelse(is.null(down), yes = "NULL", no = down), "\\.R")) %>%
    str_replace(pattern = "\\.R", replacement = str_c("_10llam", 10*as.numeric(llam_temp), "\\.R")) %>%
    str_replace(pattern = "\\.R", replacement = str_c("_PC", PC_temp, "\\.R")) 
  
  folder_file_name_R = paste(folder_name, file_name_R, sep = "/")
  
  write(Rcodetext_child,
        file = folder_file_name_R)
  
  # Change the sh file
  log_name_R = str_replace(file_name_R, pattern = "\\.R", replacement = "\\.log")
  
  Rscript_line = 
    str_c("Rscript ", file_name_R, 
          " > ", log_name_R, 
          " 2>&1")
  
  jobName_line = 
    str_replace(shcodetext[where_jobName], 
                pattern = "job-name=[[:graph:]]*", 
                replacement = str_c("job-name=", file_name_R)) %>%
    str_replace(pattern = "\\.R", replacement = "")
  
  shcodetext_child[where_Rscript] = Rscript_line
  shcodetext_child[where_jobName] = jobName_line
  
  file_name_sh = 
    str_replace(file_name_R, pattern = "\\.R", replacement = "\\.sh")
  
  folder_file_name_sh = paste(folder_name, file_name_sh, sep = "/")
  
  write(shcodetext_child, 
        file = folder_file_name_sh) 
}
}




